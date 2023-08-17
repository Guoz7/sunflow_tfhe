#include <stdlib.h>
#include <stdio.h>


char* myStrncpy(char* dest, const char* src, size_t n) {
    size_t i;
    for (i = 0; i < n && src[i] != '\0'; ++i) {
        dest[i] = src[i];
    }
    for (; i < n; ++i) {
        dest[i] = '\0';
    }
    return dest;
}


// 简化的内存块结构
typedef struct Block {
    size_t size;
    struct Block* next;
} Block;

// 简化的变量结构
typedef struct Variable {
    char name[32];
    size_t size;
    void* value;
} Variable;

// 内存池起始地址
static void* memory_start = NULL;
// 空闲内存链表
static Block* free_list = NULL;

// 初始化内存池
void initMemoryPool(void* start, size_t size) {
    memory_start = start;
    free_list = (Block*)start;
    free_list->size = size - sizeof(Block);
    free_list->next = NULL;
}

// 分配内存
void* myMalloc(size_t size) {
    if (size == 0) {
        return NULL;
    }

    Block* prev = NULL;
    Block* current = free_list;

    while (current != NULL) {
        if (current->size >= size) {
            if (current->size > size + sizeof(Block)) {
                Block* next = (Block*)((char*)current + size + sizeof(Block));
                next->size = current->size - size - sizeof(Block);
                next->next = current->next;
                current->size = size;
                current->next = next;
            }

            if (prev == NULL) {
                free_list = current->next;
            } else {
                prev->next = current->next;
            }
            printf("the address of current is %p\n", current);
            return (char*)current + sizeof(Block);
        }
        
        prev = current;
        current = current->next;
    }

    return NULL; // 没有足够的空闲内存
}

// 创建一个变量
Variable* createVariable(const char* name, size_t size) {
    Variable* var = (Variable*)myMalloc(sizeof(Variable));
    if (var != NULL) {
        var->size = size;
        var->value = myMalloc(size);
        if (var->value != NULL) {
            myStrncpy(var->name, name, sizeof(var->name) - 1);
            var->name[sizeof(var->name) - 1] = '\0';
            return var;
        } else {
            // 内存分配失败，释放之前分配的变量结构
            myFree(var);
        }
    }
    return NULL;
}

// 释放内存
void myFree(void* ptr) {
    if (ptr != NULL) {
        Block* block = (Block*)((char*)ptr - sizeof(Block));
        block->next = free_list;
        free_list = block;
    }
}

int main() {
    // 为内存池分配一块内存区域
    char memory[4096000];
    initMemoryPool(memory, sizeof(memory));

    // 创建一个变量
    Variable* var = createVariable("myVar", sizeof(int));
    if (var != NULL) {
        int value = 42;
        //memcpy(var->value, &value, sizeof(int));
        printf("the address of var->value is %p\n", &var->value);
        printf("%s = %d\n", var->name, *(int*)var->value);
    }

    // 释放变量和内存池
    myFree(var);
    
    return 0;
}
